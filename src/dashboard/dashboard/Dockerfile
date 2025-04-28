# Use an official Node.js image as the base image for building
FROM node:18-alpine AS build

# Set the working directory inside the container
WORKDIR /app

# # Copy package.json and package-lock.json before copying everything else for better caching
# COPY package.json package-lock.json ./

# # Install dependencies
# RUN npm install

# Copy the rest of the project files
COPY . .

# COPY ./taxonium_component .

# COPY ./taxonium_data_handling .

# WORKDIR /app/dashboard
# # # Build the Vite React app

# RUN yarn

# RUN yarn run build

# Use Nginx to serve the React Vite build
FROM nginx:latest AS production

# Remove the default Nginx configuration
RUN rm /etc/nginx/conf.d/default.conf

# Copy the custom Nginx configuration file
COPY ./nginx-react.conf /etc/nginx/conf.d/default.conf

# Copy the built Vite React app from the build stage
COPY --from=build /app/dist /usr/share/nginx/html

# Expose port 80 for serving the frontend
EXPOSE 80

# Start Nginx
CMD ["nginx", "-g", "daemon off;"]
